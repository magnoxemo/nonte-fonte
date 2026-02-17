file(GLOB src *.cpp *.C)

add_library(note_fonte_lib SHARED ${src})
target_link_libraries(note_fonte_lib OpenMP::OpenMP_CXX)