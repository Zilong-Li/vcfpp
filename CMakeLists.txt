cmake_minimum_required(VERSION 3.14)

project(
    vcfpp
    VERSION 0.7.4
    LANGUAGES CXX
)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

option(VCFPP_USE_PKG_CONFIG "Use pkg-config to locate htslib" ON)

add_library(vcfpp INTERFACE)
add_library(vcfpp::vcfpp ALIAS vcfpp)

target_compile_features(vcfpp INTERFACE cxx_std_11)
target_include_directories(
    vcfpp
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if(VCFPP_USE_PKG_CONFIG)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(HTSLIB QUIET IMPORTED_TARGET htslib)
    endif()
endif()

if(TARGET PkgConfig::HTSLIB)
    target_link_libraries(vcfpp INTERFACE PkgConfig::HTSLIB)
else()
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/hts.h)
    find_library(HTSLIB_LIBRARY NAMES hts)

    if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
        add_library(HTSLIB::hts UNKNOWN IMPORTED)
        set_target_properties(
            HTSLIB::hts
            PROPERTIES
                IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}"
        )
        target_link_libraries(vcfpp INTERFACE HTSLIB::hts)
    else()
        message(
            WARNING
            "htslib was not found. Consumers of vcfpp::vcfpp must make htslib headers and libhts available."
        )
        target_link_libraries(vcfpp INTERFACE hts)
    endif()
endif()

install(FILES vcfpp.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(
    TARGETS vcfpp
    EXPORT vcfppTargets
)

install(
    EXPORT vcfppTargets
    NAMESPACE vcfpp::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/vcfpp"
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/vcfppConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/vcfppConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/vcfpp"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/vcfppConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/vcfppConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/vcfppConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/vcfpp"
)
