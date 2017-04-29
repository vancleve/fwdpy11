//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
#include <fwdpy11/types.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <queue>
#include <cmath>
#include <stdexcept>
#include <fwdpp/diploid.hh>
#include <fwdpp/experimental/sample_diploid.hpp>
#include <fwdpp/experimental/sample_diploid_mloc.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpy11/samplers.hpp>
#include <fwdpy11/fitness/fitness.hpp>
#include <fwdpy11/rules/qtrait.hpp>
#include <fwdpy11/sim_functions.hpp>

namespace py = pybind11;

// Generation time, optimum, VS, sigE
using env = std::tuple<KTfwd::uint_t, double, double, double>;

static const std::size_t GEN = 0;
static const std::size_t OPTIMUM = 1;
static const std::size_t VS = 2;
static const std::size_t SIGE = 3;

// Evolve the population for some amount of time with mutation and
// recombination
void
evolve_singlepop_regions_qtrait_cpp(
    const fwdpy11::GSLrng_t &rng, fwdpy11::singlepop_t &pop,
    py::array_t<std::uint32_t> popsizes, const double mu_neutral,
    const double mu_selected, const double recrate,
    const KTfwd::extensions::discrete_mut_model &mmodel,
    const KTfwd::extensions::discrete_rec_model &rmodel,
    fwdpy11::single_locus_fitness &fitness,
    fwdpy11::singlepop_temporal_sampler recorder, const double selfing_rate,
    std::function<double(double)> trait_to_fitness,
    py::object trait_to_fitness_updater,
    std::function<double(const double g, const fwdpy11::diploid_t &,
                         const fwdpy11::diploid_t &)>
        noise,
    py::object noise_updater)
{
    bool updater_exists = false;
    py::function updater;
    if (trait_to_fitness_updater != py::none())
        {
            updater = py::function(trait_to_fitness_updater);
            updater_exists = true;
        }
    bool noise_updater_exists = false;
    py::function noise_updater_fxn;
    if (noise_updater != py::none())
        {
            noise_updater_fxn = noise_updater;
            noise_updater_exists = true;
        }
    const auto generations = popsizes.size();
    if (!generations)
        throw std::runtime_error("empty list of population sizes");
    if (mu_neutral < 0.)
        {
            throw std::runtime_error("negative neutral mutation rate: "
                                     + std::to_string(mu_neutral));
        }
    if (mu_selected < 0.)
        {
            throw std::runtime_error("negative selected mutation rate: "
                                     + std::to_string(mu_selected));
        }
    if (recrate < 0.)
        {
            throw std::runtime_error("negative recombination rate: "
                                     + std::to_string(recrate));
        }
    // if (environmental_changes.empty())
    //    {
    //        throw std::runtime_error("empty list of environmental changes.");
    //    }
    // for (auto&& ei : environmental_changes)
    //    {
    //        if (std::get<VS>(ei) < 0.)
    //            {
    //                throw std::runtime_error("VS must be >= 0.");
    //            }
    //        if (std::get<SIGE>(ei) < 0.)
    //            {
    //                throw std::runtime_error(
    //                    "environmental noise must be >= 0.");
    //            }
    //    }
    const auto fitness_callback = fitness.callback();
    pop.mutations.reserve(std::ceil(
        std::log(2 * pop.N)
        * (4. * double(pop.N) * (mu_neutral + mu_selected)
           + 0.667 * (4. * double(pop.N) * (mu_neutral + mu_selected)))));
    const auto recmap = KTfwd::extensions::bind_drm(
        rmodel, pop.gametes, pop.mutations, rng.get(), recrate);
    const auto mmodels = KTfwd::extensions::bind_dmm(
        mmodel, pop.mutations, pop.mut_lookup, rng.get(), mu_neutral,
        mu_selected, &pop.generation);
    ++pop.generation;
    auto rules = fwdpy11::qtrait::qtrait_model_rules(trait_to_fitness, noise);
    // generate FIFO queue of env changes
    // std::queue<env> env_q(std::deque<env>(environmental_changes.begin(),
    //                                      environmental_changes.end()));
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            // if (!env_q.empty())
            //    {
            //        auto next_env = env_q.front();
            //        if (pop.generation >= std::get<GEN>(next_env))
            //            {
            //                rules.optimum = std::get<OPTIMUM>(next_env);
            //                rules.VS = std::get<VS>(next_env);
            //                rules.sigE = std::get<SIGE>(next_env);
            //                env_q.pop();
            //            }
            //    }
            fitness.update(pop);
            const auto N_next = popsizes.at(generation);
            double wbar = KTfwd::experimental::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, pop.N, N_next, mu_neutral + mu_selected, mmodels,
                recmap, fitness_callback, pop.neutral, pop.selected,
                selfing_rate, rules, KTfwd::remove_neutral());
            pop.N = N_next;
            fwdpy11::update_mutations_n(
                pop.mutations, pop.fixations, pop.fixation_times,
                pop.mut_lookup, pop.mcounts, pop.generation, 2 * pop.N);
            recorder(pop);
            if (updater_exists)
                {
                    updater(pop.generation);
                }
            if (noise_updater_exists)
                {
                    noise_updater_fxn(pop.generation);
                }
        }
    --pop.generation;
}
   
//evolve_qtrait_mloc_regions_cpp(rng,pop,
//        popsizes,
//        mu_neutral,
//        mu_selected,
//        recrate,
//        mm,
//        rm,
//        interlocus_rec,
//        genetic_value_model,
//        recorder,selfing_rate,
//        multilocus_trait_model,
//        trait_to_fitness,
//        updater,noise,noise_updater)


void
evolve_qtrait_mloc_regions_cpp(
    const fwdpy11::GSLrng_t &rng, fwdpy11::multilocus_t &pop,
    py::array_t<std::uint32_t> popsizes,
    const std::vector<double> &neutral_mutation_rates,
    const std::vector<double> &selected_mutation_rates,
    const std::vector<double> &recrates,
    const std::vector<KTfwd::extensions::discrete_mut_model> &mmodels,
    const std::vector<KTfwd::extensions::discrete_rec_model> &rmodels,
    const std::vector<std::function<unsigned(void)>> & interlocus_rec,
    fwdpy11::multilocus_genetic_value &multilocus_gvalue,
    fwdpy11::multilocus_temporal_sampler recorder, const double selfing_rate,
    std::function<double(const py::array_t<double>)> aggregator,
    std::function<double(double)> trait_to_fitness,
    py::object trait_to_fitness_updater,
    std::function<double(const double g, const fwdpy11::multilocus_diploid_t &,
                         const fwdpy11::multilocus_diploid_t &)>
        noise,
    py::object noise_updater)
{
    bool updater_exists = false;
    py::function updater;
    if (trait_to_fitness_updater != py::none())
        {
            updater = py::function(trait_to_fitness_updater);
            updater_exists = true;
        }
    bool noise_updater_exists = false;
    py::function noise_updater_fxn;
    if (noise_updater != py::none())
        {
            noise_updater_fxn = noise_updater;
            noise_updater_exists = true;
        }
    const auto generations = popsizes.size();
    if (!generations)
        throw std::runtime_error("empty list of population sizes");
    auto bound_mmodels = KTfwd::extensions::bind_vec_dmm(
        mmodels, pop.mutations, pop.mut_lookup, rng.get(),
        neutral_mutation_rates, selected_mutation_rates, &pop.generation);
    auto bound_intralocus_rec = KTfwd::extensions::bind_vec_drm(
        rmodels, pop.gametes, pop.mutations, rng.get(), recrates);
    std::vector<double> total_mut_rates(neutral_mutation_rates);
    std::transform(total_mut_rates.cbegin(), total_mut_rates.cend(),
                   selected_mutation_rates.cbegin(), total_mut_rates.begin(),
                   std::plus<double>());

    fwdpy11::qtrait::qtrait_mloc_rules rules(aggregator, trait_to_fitness,
                                             noise);

    ++pop.generation;
    for (unsigned i = 0; i < generations; ++i, ++pop.generation)
        {
            auto N_next = popsizes.at(i);
            auto wbar = KTfwd::experimental::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, pop.N, N_next, total_mut_rates.data(),
                bound_mmodels, bound_intralocus_rec, interlocus_rec,
                multilocus_gvalue, pop.neutral, pop.selected, selfing_rate,
                rules, KTfwd::remove_neutral());
            pop.N = N_next;
            fwdpy11::update_mutations_n(
                pop.mutations, pop.fixations, pop.fixation_times,
                pop.mut_lookup, pop.mcounts, pop.generation, 2 * pop.N);
            recorder(pop);
            if (updater_exists)
                {
                    updater(pop.generation);
                }
            if (noise_updater_exists)
                {
                    noise_updater_fxn(pop.generation);
                }
            --pop.generation;
        }
}

PYBIND11_PLUGIN(wfevolve_qtrait)
{
    py::module m("wfevolve_qtrait", "example extending");

    m.def("evolve_singlepop_regions_qtrait_cpp",
          &evolve_singlepop_regions_qtrait_cpp);

    m.def("evolve_qtrait_mloc_regions_cpp",
          &evolve_qtrait_mloc_regions_cpp);
    return m.ptr();
}
