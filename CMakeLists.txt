cmake_minimum_required(VERSION 3.20)
project(Enterprise C)

set(CMAKE_C_STANDARD 11)

add_executable(Enterprise
        main.c
        orbit_types.h
        astrodynamics.h
        astrodynamics.c
        astrodynamics.c
        propagator.h
        propagator.c
        telemetry.h
        telemetry.c
        mission_analysis.h
        mission_analysis.c
)