CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

SET(EXECUTABLE_NAME morph_grid)

SET(SRC
	bilateral_filter_future.cpp
	bilateral_filter_future.hpp
	canny_future.cpp
	canny_future.hpp
	edge.hpp
	ez_grid.cpp
	ez_grid.hpp
	future_base.hpp
	morph_grid.cpp
	morph_grid_future.cpp
	morph_grid_future.hpp
	${PROJECT_SOURCE_DIR}/config.h
)

ADD_EXECUTABLE(${EXECUTABLE_NAME} ${SRC})

TARGET_LINK_LIBRARIES(${EXECUTABLE_NAME} ${LIBS})