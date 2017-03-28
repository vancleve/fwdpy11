def makeMutationRegions(neutral,selected):
    """
    Convert user input into :class:`~fwdpy11.fwdpp_extensions.MutationRegions`

    :param neutral: A list of :class:`fwdpy11.regions.Region` objects.
    :param selected: A list of :class:`fwdpy11.regions.Sregion` objects.

    :rtype: :class:`fwdpy11.fwdpp_extensions.MutationRegions`

    .. note:: Used by various "evolve" functions.  Users probably won't need to call this.
    
        >>> import fwdpy11 as fp11
        >>> nregions = [fp11.Region(0,0.5,1),fp11.Region(1,1.5,1)]
        >>> sregions = [fp11.ExpS(0,0.5,1,1),fp11.GaussianS(1,1.5,1,0.25)]
        >>> mr = fp11.makeMutationRegions(nregions,sregions)
        >>> type(mr)
        <class 'fwdpy11.fwdpp_extensions.MutationRegions'>

    One or both lists may be empty:

        >>> mr = fp11.makeMutationRegions([],[])

    Neither list may be None:

        >>> mr = fp11.makeMutationRegions([],None)
        Traceback (most recent call last):
            ...
        TypeError: 'NoneType' object is not iterable
    """
    nbeg = [i.b for i in neutral]
    nend = [i.e for i in neutral]
    nweights = [i.w for i in neutral]
    sbeg = [i.b for i in selected]
    send = [i.e for i in selected]
    sweights = [i.w for i in selected]
    sh = [i.callback() for i in selected]
    from .fwdpp_extensions import MutationRegions
    return MutationRegions(nbeg,nend,nweights,sbeg,send,sweights,sh)

def makeRecombinationRegions(regions):
    """
    Convert user input into :class:`~fwdpy11.fwdpp_extensions.RecombinationRegions`

    :param neutral: A list of :class:`fwdpy11.regions.Region` objects.
    :param selected: A list of :class:`fwdpy11.regions.Sregion` objects.

    :rtype: :class:`fwdpy11.fwdpp_extensions.MutationRegions`

    .. note:: Used by various "evolve" functions.  Users probably won't need to call this.
    
        >>> import fwdpy11 as fp11
        >>> rregions = [fp11.Region(0,0.5,1),fp11.Region(1,1.5,1)]
        >>> rr = fp11.makeRecombinationRegions(rregions)
        >>> type(rr)
        <class 'fwdpy11.fwdpp_extensions.RecombinationRegions'>

    """
    beg = [i.b for i in regions]
    end = [i.e for i in regions]
    weights = [i.w for i in regions]
    from .fwdpp_extensions import RecombinationRegions
    return RecombinationRegions(beg,end,weights)
