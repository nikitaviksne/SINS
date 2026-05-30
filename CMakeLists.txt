cmake_minimum_required(VERSION 3.11)

project("Inertial_algorithm" LANGUAGES CXX)

set(SRCS
    mathematics.cpp
    matrix.cpp
    )

set(SOURCE_LIB
    include/
	/home/nikita/Документы/CLibraries/alglib-cpp-v4.04.0/src/
    )

add_library(lclib STATIC ${SRCS})
target_include_directories(lclib PUBLIC ${SOURCE_LIB})

# Библиотека с алгоритмами БИНС
set(Headers_SINS 
	SINS_Algorithms/
	include/
	/home/nikita/Документы/CLibraries/alglib-cpp-v4.04.0/src/		
	)
#message(STATUS Headers_SINS, ${Headers_SINS})
set(Source_SINS 
		SINS_Algorithms/alignment.cpp
		SINS_Algorithms/SolveOrient.cpp
		SINS_Algorithms/SolveNav.cpp
		mathematics.cpp
		matrix.cpp
	)

#message(STATUS Source_SINS, ${Source_SINS})

add_library(ins STATIC ${Source_SINS}) #Была STATIC

target_include_directories(ins PUBLIC ${Headers_SINS})

#для Калмана требуется библиотека alglib, она у меня уже собрана, надо ее найти
# find_library(alglib NAMES libalglib.a PATHS "Kalman/") #не сработало
# set(Headers_alglib
# 		../../CLibraries/alglib-cpp-v4.04.0/src/
# 	)

# set(Source_alglib
# 		../../CLibraries/alglib-cpp-v4.04.0/src/alglibmisc.cpp
# 		../../CLibraries/alglib-cpp-v4.04.0/src/alglibinternal.cpp
# 		../../CLibraries/alglib-cpp-v4.04.0/src/linalg.cpp
# 		../../CLibraries/alglib-cpp-v4.04.0/src/ap.cpp	
# 	)
# add_library(alglib STATIC ${Source_alglib})
# target_include_directories(alglib PUBLIC ${Headers_alglib})

#добавляем библиотеку alglib, как статическую
set(SOURCE_ALGLIB
		/home/nikita/Документы/CLibraries/alglib-cpp-v4.04.0/src/
    )

add_library(alglib STATIC 
		/home/nikita/Документы/CLibraries/alglib-cpp-v4.04.0/src/linalg.cpp
    	/home/nikita/Документы/CLibraries/alglib-cpp-v4.04.0/src/ap.cpp
    	/home/nikita/Документы/CLibraries/alglib-cpp-v4.04.0/src/alglibinternal.cpp
    	/home/nikita/Документы/CLibraries/alglib-cpp-v4.04.0/src/alglibmisc.cpp
	)
target_include_directories(alglib PUBLIC ${SOURCE_ALGLIB})# указываем заголовочные файлы


# Библиотека фильтра Калмана
set(Headers_Kalman
		Kalman/
		include/
		/home/nikita/Документы/CLibraries/alglib-cpp-v4.04.0/src/		
	)
set(Source_Kalman 
		Kalman/AdaptiveRKalman.cpp
		#Kalman/AdaptiveQRKalman.cpp
		Kalman/UsualKalman.cpp
		Kalman/mathematics_kalman.cpp	
	)
#message(STATUS Source_Kalman, ${Source_Kalman})
add_library(kalman STATIC ${Source_Kalman})
target_include_directories(kalman PUBLIC ${Headers_Kalman})

# Библиотека кватернионов
set(Headers_Quaternions
		include/
	)
set(Source_Quaternions
		quaternions.cpp
	)
#message(STATUS Source_Kalman, ${Source_Kalman})
add_library(quaternions STATIC ${Source_Quaternions})
target_include_directories(quaternions PUBLIC ${Headers_Quaternions})

# Tell CMake to run moc when necessary:
#set(CMAKE_AUTOMOC ON)
#set(CMAKE_AUTOUIC ON)
# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
#set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Widgets finds its own dependencies.
#find_package(Qt5Widgets REQUIRED) #Выдает ошибку о требовании GLIBCXX_3.4.32 если компилировать вместе с alblib.

add_executable(main main.cpp)
target_link_directories(main PUBLIC 
						include/)
#qt5_use_modules(main Core) ##Выдает ошибку о требовании GLIBCXX_3.4.32 если компилировать вместе с alblib.

target_link_libraries(main 
	lclib 
	ins 
	kalman 
	alglib
	quaternions
)
