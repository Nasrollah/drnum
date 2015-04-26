SET(testBlockObjects_CC_SOURCES
compressibleeulerlsobc.cpp
compressibleeulerlsobc.h
compressiblesimpleswalllsobc.cpp
compressiblesimpleswalllsobc.h
compressibleswalllsobc.cpp
compressibleswalllsobc.h
geoblockobjecttest001.h
geoblockobjecttest002.h
geoblockobjecttest003.h
geoblockobjecttest004.h
geoblockobjecttest005.h
geoblockobjecttest006.h
geolevelsettest001.h
geolevelsettest002.h
geolevelsettest007.h
lsbc_cell_t.h
lsbc_list_t.h
lsobject.h
lsobject.cpp
main.cpp
main.h
main.cu)

SET(DRNUM_USED_LIBS drnumlib shmlib)

IF(USE_GPU)
    SET(testBlockObjects_CU_SOURCES main.cu)
    cuda_add_library(testBlockObjects_cuda ${testBlockObjects_CU_SOURCES} SHARED) # OPTIONS -G)
    SET(DRNUM_USED_LIBS testBlockObjects_cuda ${DRNUM_USED_LIBS})
ENDIF(USE_GPU)


ADD_EXECUTABLE(testBlockObjects ${testBlockObjects_CC_SOURCES})
ADD_DEPENDENCIES(testBlockObjects ${DRNUM_USED_LIBS})
TARGET_LINK_LIBRARIES(testBlockObjects ${DRNUM_USED_LIBS} ${QT_LIBRARIES} ${VTK_LIBRARIES} ${MPI_LIBRARIES} ${OPENMP_LIBS})
#QT5_USE_MODULES(drnumJetDemo Core)

SET_TARGET_PROPERTIES(testBlockObjects
    PROPERTIES
    LINKER_LANGUAGE CXX
    PREFIX "")

#SET_TARGET_PROPERTIES(testBlockObjects
#    PROPERTIES
#    VERSION ${DRNUM_VERSION})

INSTALL(TARGETS testBlockObjects RUNTIME DESTINATION bin)
#IF(USE_GPU)
#    INSTALL(TARGETS drnumJetDemo_cuda LIBRARY DESTINATION lib)
#ENDIF(USE_GPU)

