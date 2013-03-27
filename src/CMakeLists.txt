file(GLOB Tyche_SOURCES "*.cpp")
include_directories(${Tyche_INCLUDE_DIRECTORIES})

add_library(Tyche SHARED ${Moirai_SOURCES})

