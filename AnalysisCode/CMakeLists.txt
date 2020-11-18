include (PKG)
PKG(AnalysisCode
    DEPENDS
        PMTSim
        SimUtil
        DetSimAlg

        GenTools

        SimEventV2
        EvtNavigator

        boost_filesystem boost_system
        Geom
)
