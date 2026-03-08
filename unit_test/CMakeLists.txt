file(GLOB tests *.cxx *.cpp)


foreach(test ${tests})
    get_filename_component(test_name ${test} NAME_WE)
    add_executable(${test_name} ${test})
    target_link_libraries(${test_name}
            PRIVATE
            note_fonte_lib
            GTest::GTest
            GTest::Main

    )
endforeach()