/**
 *  @file    gitversion.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2018-06-30
 *
 *  @brief This file together with a makefile entry lets us access current git revision
 *
 *  The makefile entry should look something like this:
 * 
 *  src/gitversion.cc: ../.git/HEAD ../.git/index $(filter-out src/gitversion.cc, $(CXX_FILES))                                                                 
 *      echo "const char* gitversion = \"$(shell git rev-parse HEAD)$(shell git diff --quiet || echo '-dirty')\";" > $@
 * 
 */

#pragma once

extern const char* gitversion;
