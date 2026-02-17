add_executable(simple_pdf simple_pdf.cpp)
add_executable(fet_one_dimension one_dimensional_fet.C)

target_link_libraries(simple_pdf
        PRIVATE
        note_fonte_lib
)


target_link_libraries(fet_one_dimension
        PRIVATE
        note_fonte_lib
)

