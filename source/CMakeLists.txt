add_executable(naivecounting naive_counting.cpp)
target_link_libraries (naivecounting seqan3::seqan3)

add_executable(flajoletmartin flajolet_martin.cpp)
target_link_libraries (flajoletmartin seqan3::seqan3)

add_executable(hyperloglog hyperloglog.cpp)
target_link_libraries (hyperloglog seqan3::seqan3)

add_executable(similarity similarity.cpp)
target_link_libraries (similarity seqan3::seqan3)