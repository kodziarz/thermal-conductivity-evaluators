file(READ "${INPUT}" KERNEL_SOURCE)

file(WRITE "${OUTPUT}" "#pragma once\n")
file(APPEND "${OUTPUT}" "#include <string_view>\n")
file(APPEND "${OUTPUT}" "namespace kernels {\n")
file(APPEND "${OUTPUT}" "inline constexpr std::string_view ${VARIABLE} = R\"CLC(")
file(APPEND "${OUTPUT}" "${KERNEL_SOURCE}")
file(APPEND "${OUTPUT}" ")CLC\";\n")
file(APPEND "${OUTPUT}" "}\n")
