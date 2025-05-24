cmake_minimum_required(VERSION 3.16)

project(newton-fb-QT VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)

# أضف ملفات muParser هنا
set(MUPARSER_SRC
    muparser/muParser.cpp
    muparser/muParserBase.cpp
    muparser/muParserBytecode.cpp
    muparser/muParserCallback.cpp
    muparser/muParserDLL.cpp
    muparser/muParserError.cpp
    muparser/muParserInt.cpp
    muparser/muParserTokenReader.cpp
)

set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
        mainwindow.ui
        ${MUPARSER_SRC} # أضف ملفات muParser لمصادر المشروع
)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(newton-fb-QT
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
    )
else()
    if(ANDROID)
        add_library(newton-fb-QT SHARED
            ${PROJECT_SOURCES}
        )
    else()
        add_executable(newton-fb-QT
            ${PROJECT_SOURCES}
        )
    endif()
endif()

# تضمين ملفات muParser
target_include_directories(newton-fb-QT PRIVATE muparser)

# تحديد تعريف ثابت لمكتبة muParser
target_compile_definitions(newton-fb-QT PRIVATE MUPARSER_STATIC)

# ربط Qt
target_link_libraries(newton-fb-QT PRIVATE Qt${QT_VERSION_MAJOR}::Widgets)



if(${QT_VERSION} VERSION_LESS 6.1.0)
  set(BUNDLE_ID_OPTION MACOSX_BUNDLE_GUI_IDENTIFIER com.example.newton-fb-QT)
endif()
set_target_properties(newton-fb-QT PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

include(GNUInstallDirs)
# install(newton-fb-QT
#     BUNDLE DESTINATION .
#     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
# )

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(newton-fb-QT)
endif()
