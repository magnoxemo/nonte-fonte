file(GLOB EXAMPLES *.cxx *.cpp)

# adding test executables for every case
foreach(EXAMPLE ${EXAMPLES})
    get_filename_component(EX_NAME ${EXAMPLE} NAME_WE)
    add_executable(${EX_NAME} ${EXAMPLE})
    target_link_libraries(${EX_NAME}
            PRIVATE
            note_fonte_lib
    )
endforeach()