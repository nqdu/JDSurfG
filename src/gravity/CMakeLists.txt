file(GLOB SOURCE "*.cpp")
#message("myfile = ${SOURCE}")
#file(REMOVE ${SOURCE} main.cpp )
list(FILTER SOURCE EXCLUDE REGEX "main.cpp")

#message("myfile = ${SOURCE}")

add_library(gravity ${SOURCE})