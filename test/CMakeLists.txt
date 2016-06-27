add_executable(hbsreadtest HBSReadTest.cpp)
add_executable(bbtest BoundingBoxTest.cpp)
add_executable(bernsteintest BernsteinPolynomialTest.cpp)

target_link_libraries(hbsreadtest confnurbs)
target_link_libraries(bbtest confnurbs )
target_link_libraries(bernsteintest confnurbs)

add_test(hbsreadtest hbsreadtest ${confnurbs_SOURCE_DIR}/input/sphere.dat)
add_test(boundingboxtest bbtest ${confnurbs_SOURCE_DIR}/input/sphere.dat)
add_test(bernsteintest bernsteintest)
