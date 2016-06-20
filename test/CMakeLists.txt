add_executable(hbsreadtest HBSReadTest.cpp)
add_executable(bbtest BoundingBoxTest.cpp)

target_link_libraries(hbsreadtest confnurbs)
target_link_libraries(bbtest confnurbs )

add_test(test1 hbsreadtest ${confnurbs_SOURCE_DIR}/input/sphere.dat)
add_test(test2 bbtest ${confnurbs_SOURCE_DIR}/input/sphere.dat)
