add_executable(simple_pdf simple_pdf.cpp)

target_link_libraries(simple_pdf
        PRIVATE
        note_fonte_lib
)

