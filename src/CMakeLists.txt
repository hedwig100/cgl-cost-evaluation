
include_directories(${PROJECT_SOURCE_DIR}/lib)

## drawSSIG
add_executable(drawSSIG
    drawSSIG.cpp
)

## cgl
add_executable(cgl
    cgl.cpp
)

## measure_op_cost
add_executable(measure_op_cost
    measure_op_cost.cpp
)

## modular_cgl
add_executable(modular_cgl
    modular_cgl.cpp
)

## mpCGL
add_executable(hn23
    hn23.cpp
)

## mpCGL
add_executable(mpCGL
    mpCGL.cpp
)

## mycgl
add_executable(mycgl
    mycgl.cpp
)

## count-cgl
add_executable(count-cgl
    count-cgl.cpp
)

## dpb_cgl
add_executable(dpb_cgl
    dpb_cgl.cpp
)

## dpb_collision
add_executable(dpb_collision
    dpb_collision.cpp
)

## find_collision
add_executable(find_collision
    find_collision.cpp
)

## radical_cgl
add_executable(radical_cgl
    radical_cgl.cpp
)

## safe_dpb_cgl
add_executable(safe_dpb_cgl
    safe_dpb_cgl.cpp
)

## time-cgl
add_executable(time-cgl
    time-cgl.cpp
)