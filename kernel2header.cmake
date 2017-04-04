#http://www.scottenglert.com/free-stuff/adding-opencl-kernel-to-gpu-deformer-without-kernel-source-file/

set(include_file "${CMAKE_ARGV3}/${CMAKE_ARGV4}")
set(kernel_file "${CMAKE_ARGV3}/${CMAKE_ARGV5}")
set(header_file "${CMAKE_ARGV3}/${CMAKE_ARGV6}")
 
file(READ "${include_file}" includeCode)
string(REGEX REPLACE ";" "\\\\;" includeCode "${includeCode}")
string(REGEX REPLACE "\"" "\\\\\"" includeCode "${includeCode}")
#added to handle blank lines when iterating
string(REGEX REPLACE "\n" "E;" includeCode "${includeCode}")
 
# Open header file to write to
file(WRITE "${header_file}" "std::string kernelSourceStr = (\"\\n\" \\\n")
#file(APPEND "${header_file}" "// #include boundary.h")

foreach(clLine ${includeCode})
    # Get rid of the trailing 'E'
    string(REGEX REPLACE "^(.*)E$" "\\1" line "${clLine}")
    file(APPEND "${header_file}" "    \"${line}\\n\" \\\n")
endforeach()

file(APPEND "${header_file}" "\"\\n\" \\\n")

file(READ "${kernel_file}" kernelCode)
string(REGEX REPLACE ";" "\\\\;" kernelCode "${kernelCode}")
string(REGEX REPLACE "\"" "\\\\\"" kernelCode "${kernelCode}")
#added to handle blank lines when iterating
string(REGEX REPLACE "\n" "E;" kernelCode "${kernelCode}")

foreach(clLine ${kernelCode})
    # Get rid of the trailing 'E'
    string(REGEX REPLACE "^(.*)E$" "\\1" line "${clLine}")
    file(APPEND "${header_file}" "    \"${line}\\n\" \\\n")
endforeach()
 
# write null terminated line
file(APPEND "${header_file}" "\"\\n\");\n")
