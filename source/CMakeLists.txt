add_executable(count count.cpp)
target_link_libraries (count seqan3::seqan3)

add_executable(similarity similarity.cpp)
target_link_libraries (similarity seqan3::seqan3)