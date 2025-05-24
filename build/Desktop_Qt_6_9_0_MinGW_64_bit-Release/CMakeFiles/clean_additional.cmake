# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "CMakeFiles\\newton-fb-QT_autogen.dir\\AutogenUsed.txt"
  "CMakeFiles\\newton-fb-QT_autogen.dir\\ParseCache.txt"
  "newton-fb-QT_autogen"
  )
endif()
