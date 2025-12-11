add_library(Matrix INTERFACE)
target_compile_features(Matrix INTERFACE cxx_std_23)
target_include_directories(Matrix INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)