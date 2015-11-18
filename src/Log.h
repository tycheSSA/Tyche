/*
 * MyAssert.h
 * 
 * Copyright 2012 Martin Robinson
 *
 * This file is part of RD_3D.
 *
 * RD_3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RD_3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with RD_3D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 28 Nov 2012
 *      Author: robinsonm
 */

#ifndef LOG_H_
#define LOG_H_

#ifdef _WIN32
#include <winbase.h>
#else
#include <signal.h>
#endif

#ifdef _WIN32
#define DEBUG_BREAK DebugBreak()
#else
#define DEBUG_BREAK raise(SIGTRAP)
#endif
 

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            DEBUG_BREAK; \
       } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

#define CHECK(condition, message) \
		if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            DEBUG_BREAK; \
        }

#define ERROR(message) \
            std::cerr << "Error at " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            DEBUG_BREAK;


#ifndef LOG_LEVEL
#	ifdef NDEBUG
#		define LOG_LEVEL 1
#	else
#		define LOG_LEVEL 2
#	endif
#endif


#define LOG(level, message) \
    if (level <= LOG_LEVEL) { \
    	std::cout << message << std::endl; \
    }



#endif /* LOG_H_ */
