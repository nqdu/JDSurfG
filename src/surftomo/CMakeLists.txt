file(GLOB SOURCE "*.cpp")
list(FILTER SOURCE EXCLUDE REGEX "main.cpp")

add_library(surftomo ${SOURCE})