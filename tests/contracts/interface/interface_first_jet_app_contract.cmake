cmake_minimum_required(VERSION 3.20)
if(NOT EXISTS "${LAMP_EXE}")
    message(FATAL_ERROR "Built LAMP executable required")
endif()
set(ENV{OMP_NUM_THREADS} 2)
set(ENV{OMP_THREAD_LIMIT} 2)
set(ENV{OMP_DYNAMIC} false)
set(ENV{OPENBLAS_NUM_THREADS} 1)
set(ENV{MKL_NUM_THREADS} 1)
string(TIMESTAMP _stamp "%Y%m%d_%H%M%S")
string(RANDOM LENGTH 8 ALPHABET 0123456789abcdef _suffix)
set(_run "${RUN_ROOT}/${_stamp}_${_suffix}")
file(MAKE_DIRECTORY "${_run}")
file(WRITE "${RUN_ROOT}/last_run.txt" "${_run}")
set(_fixtures "${CMAKE_CURRENT_LIST_DIR}/../../fixtures/interface_transport")
function(write_case name fixture)
    file(MAKE_DIRECTORY "${_run}/${name}")
    file(COPY "${_fixtures}/${fixture}/system" DESTINATION "${_run}/${name}")
endfunction()
function(replace_input name file old new)
    set(path "${_run}/${name}/system/${file}.dict")
    file(READ "${path}" body)
    string(FIND "${body}" "${old}" index)
    if(index LESS 0)
        message(FATAL_ERROR "Missing mutation anchor: ${file}: ${old}")
    endif()
    string(REPLACE "${old}" "${new}" body "${body}")
    file(WRITE "${path}" "${body}")
endfunction()
function(check_input name success expected)
    execute_process(COMMAND "${LAMP_EXE}" --case "${_run}/${name}" --check-case
        RESULT_VARIABLE rc OUTPUT_FILE "${_run}/${name}/check.log" ERROR_FILE "${_run}/${name}/check-error.log")
    file(READ "${_run}/${name}/check.log" out)
    file(READ "${_run}/${name}/check-error.log" err)
    if(success)
        if(NOT rc EQUAL 0 OR NOT out MATCHES "${expected}")
            message(FATAL_ERROR "Valid input rejected: ${name}, ${rc}, ${err}")
        endif()
    elseif(rc EQUAL 0 OR NOT "${out}${err}" MATCHES "${expected}")
        message(FATAL_ERROR "Expected precise input rejection absent: ${name}, ${rc}, ${out}${err}")
    endif()
    if(out MATCHES "Building background mesh" OR EXISTS "${_run}/${name}/output" OR EXISTS "${_run}/${name}/cache")
        message(FATAL_ERROR "check-case generated mesh/output: ${name}")
    endif()
    file(WRITE "${_run}/${name}/check-exit.txt" "${rc}\n")
endfunction()
set(_valid 0)
foreach(fixture static_xy static_rz translation_xy translation_rz rotation_xy vortex_xy)
    foreach(mesh uniform amr)
        set(name "${fixture}_${mesh}")
        write_case("${name}" "${fixture}")
        replace_input("${name}" mesh "initialMesh uniform;" "initialMesh ${mesh};\nRefinement { min (.25 .25); max (.75 .75); level 1; }\nRegrid { atSteps (1 2 3); }")
        if(mesh STREQUAL amr)
            foreach(file geometry interface)
                file(READ "${_run}/${name}/system/${file}.dict" body)
                string(REGEX REPLACE "phi \"([^\"]+)\"" "phi \"exp(\\1)-1\"" body "${body}")
                file(WRITE "${_run}/${name}/system/${file}.dict" "${body}")
            endforeach()
        endif()
        check_input("${name}" TRUE "FIRST_JET_APP_INPUT_OK")
        execute_process(COMMAND "${LAMP_EXE}" --case "${_run}/${name}"
            RESULT_VARIABLE rc OUTPUT_FILE "${_run}/${name}/run.log" ERROR_FILE "${_run}/${name}/error.log")
        file(WRITE "${_run}/${name}/exit.txt" "${rc}\n")
        file(READ "${_run}/${name}/run.log" out)
        if(NOT rc EQUAL 0 OR NOT out MATCHES "FIRST_JET_APP_RESULT,(xy|rz),8,[^\n]+,3,COMPLETED")
            message(FATAL_ERROR "Actual five-file FirstJet smoke failed: ${name}, ${rc}; ${_run}")
        endif()
        string(REGEX MATCHALL "\\[Step [0-9]+\\]" accepted "${out}")
        list(LENGTH accepted count)
        string(REGEX MATCHALL "FIRST_JET_APP_REGRID," regrids "${out}")
        list(LENGTH regrids regrid_count)
        if(NOT count EQUAL 4 OR NOT regrid_count EQUAL 3 OR NOT out MATCHES "inventory=NA flux_ledger=NA")
            message(FATAL_ERROR "Missing real accepted/regrid/NA records: ${name}")
        endif()
        math(EXPR _valid "${_valid}+1")
    endforeach()
endforeach()
set(_invalid 0)
function(reject_mutation name fixture file old new reason)
    write_case("${name}" "${fixture}")
    replace_input("${name}" "${file}" "${old}" "${new}")
    check_input("${name}" FALSE "${reason}")
    math(EXPR count "${_invalid}+1")
    set(_invalid ${count} PARENT_SCOPE)
endfunction()
foreach(spec
    "unknown|interface|method hermiteGals;|method hermiteGals; unused 1;|unknown input"
    "duplicate|run|dt 0.005;|dt 0.005; dt .001;|duplicate input"
    "problem|interface|method hermiteGals;|method hermiteGals; problem staticCircleXY;|unknown input"
    "mesh_word|mesh|initialMesh uniform|initialMesh adaptive|uniform.*amr"
    "boolean|run|enabled false|enabled TRUE|true.*false"
    "method|interface|method hermiteGals|method unknown|Interface.method.*expected hermiteGals"
    "unused_schedule|run|enabled false;|enabled false; writeEverySteps 1;|unused schedule"
    "unused_constant|fields|radius .19;|radius .19; unused 1;|unused constant"
    "reserved_constant|fields|radius .19;|radius .19; t 1;|nonreserved"
    "coordinates|geometry|coordinates xy|coordinates xyz|xy.*rz"
    "bad_formula|geometry|(x-cx)*(x-cx)|undefined(x)|Geometry.Initial.phi.*valid analytic expression"
    "initial_time|geometry|(x-cx)*(x-cx)|x+t|Geometry.Initial.phi.*unknown identifier"
    "singular|geometry|(x-cx)*(x-cx)|1/x|Geometry.Initial.phi.*division by zero"
    "unused_boundary|interface|type prescribed|type none|unknown input"
    "unused_velocity|fields|type constant;|type constant; origin (0 0);|unknown input"
    "vector|fields|value (0 0)|value (0)|two-component"
    "domain|geometry|max (1 1)|max (0 1)|ordered bounds"
    "rz_domain|geometry|coordinates xy|coordinates rz|Geometry.Initial.phi.*unknown identifier"
    "junk|run|dt 0.005|dt .005junk|finite number"
    "nan|run|dt 0.005|dt nan|finite number"
    "units|run|units code|units SI|expected code"
    "fractional_window|run|endTime 0.02|endTime .0201|integer number"
    "two_clocks|run|dt 0.005;|dt 0.005; steps 4;|unknown input"
    "mesh_in_interface|interface|method hermiteGals;|method hermiteGals; resolution 8;|Interface.resolution.*unknown input"
    "time_in_interface|interface|method hermiteGals;|method hermiteGals; dt .01;|Interface.dt.*unknown input"
    "geometry_in_interface|interface|method hermiteGals;|method hermiteGals; Domain { min (0 0); max (1 1); }|Interface.Domain.*unknown input"
    "velocity_in_interface|interface|method hermiteGals;|method hermiteGals; Velocity { type constant; value (0 0); }|Interface.Velocity.*unknown input"
    "output_in_interface|interface|method hermiteGals;|method hermiteGals; writeFields true;|Interface.writeFields.*unknown input"
    "constant_in_interface|interface|method hermiteGals;|method hermiteGals; Constants { radius .19; }|Interface.Constants.*unknown input"
    "initial_in_interface|interface|method hermiteGals;|method hermiteGals; Initial { phi x; }|Interface.Initial.*unknown input"
    "mesh_in_run|run|Run { app Interface; }|Run { app Interface; resolution 8; }|Run.resolution.*unknown input"
    "method_in_mesh|mesh|resolution 8;|resolution 8; method hermiteGals;|Mesh.method.*unknown input")
    string(REGEX MATCH "^([^|]+)\\|([^|]+)\\|([^|]+)\\|([^|]+)\\|(.+)$" matched "${spec}")
    reject_mutation("${CMAKE_MATCH_1}" static_xy "${CMAKE_MATCH_2}" "${CMAKE_MATCH_3}" "${CMAKE_MATCH_4}" "${CMAKE_MATCH_5}")
endforeach()
set(_context 0)
foreach(spec "mesh|resolution 8|resolution 3|Mesh.resolution|resolution must" "mesh|resolution 8|resolution 257|Mesh.resolution|resolution must"
    "run|dt 0.005|dt 0|Time.dt|finite and positive" "run|dt 0.005|dt -1|Time.dt|finite and positive"
    "run|endTime 0.02|endTime 0|Time.endTime|finite and positive" "run|endTime 0.02|endTime 1e308|Time.endTime|integer number"
    "geometry|max (1 1)|max (0 1)|Geometry.Domain|ordered bounds")
    string(REPLACE "|" ";" parts "${spec}")
    list(GET parts 0 file)
    list(GET parts 1 old)
    list(GET parts 2 new)
    list(GET parts 3 field)
    list(GET parts 4 reason)
    reject_mutation("range_${_context}" static_xy "${file}" "${old}" "${new}" "${reason}")
    file(READ "${_run}/range_${_context}/system/${file}.dict" body)
    if(field STREQUAL "Geometry.Domain")
        set(anchor "Domain")
    else()
        set(anchor "${new}")
    endif()
    string(FIND "${body}" "${anchor}" index)
    string(SUBSTRING "${body}" 0 ${index} prefix)
    string(REGEX MATCHALL "\n" lines "${prefix}")
    list(LENGTH lines line)
    math(EXPR line "${line}+1")
    file(READ "${_run}/range_${_context}/check-error.log" err)
    string(FIND "${err}" "system/${file}.dict:${line}: ${field} = '" context)
    if(context LESS 0)
        message(FATAL_ERROR "Missing precise cross-file source context: ${err}")
    endif()
    math(EXPR _context "${_context}+1")
endforeach()
foreach(ids "3 1" "1 1" "0 2" "4")
    reject_mutation("regrid_${_invalid}" static_xy mesh "initialMesh uniform;" "initialMesh uniform; Refinement { min (.25 .25); max (.75 .75); level 1; } Regrid { atSteps (${ids}); }" "increase strictly|distinct positive|before the final")
endforeach()
foreach(spec "scale (1 1)|scale (0 1)|positive polynomial" "c01|c41|c00..c33" "-omega|x|Fields.Velocity.Coefficients.c01.x.*unknown identifier")
    string(REPLACE "|" ";" parts "${spec}")
    list(GET parts 0 old)
    list(GET parts 1 new)
    list(GET parts 2 reason)
    reject_mutation("polynomial_${_invalid}" rotation_xy fields "${old}" "${new}" "${reason}")
endforeach()
write_case(old_schema static_xy)
file(WRITE "${_run}/old_schema/system/interface.dict" "Interface { Benchmark { method hermiteGals; problem solidRotationXY; } }")
check_input(old_schema FALSE "five-file Interface input")
foreach(method levelSet conservativeGals parametricGals)
    reject_mutation("legacy_${method}" static_xy interface "method hermiteGals" "method ${method}" "Interface.method.*expected hermiteGals")
endforeach()
set(_missing 0)
foreach(file geometry mesh fields interface)
    write_case("missing_${file}" static_xy)
    file(REMOVE "${_run}/missing_${file}/system/${file}.dict")
    check_input("missing_${file}" FALSE "system/${file}.dict.*cannot open dictionary")
    math(EXPR _missing "${_missing}+1")
endforeach()
write_case(extra_file static_xy)
file(WRITE "${_run}/extra_file/system/unused.dict" "Unused { x 1; }\n")
check_input(extra_file FALSE "unused standalone Interface input")
write_case(cfl translation_xy)
replace_input(cfl run "dt 0.01" "dt 1")
replace_input(cfl run "endTime 0.040000000000000001" "endTime 2")
check_input(cfl TRUE "FIRST_JET_APP_INPUT_OK")
execute_process(COMMAND "${LAMP_EXE}" --case "${_run}/cfl" RESULT_VARIABLE rc OUTPUT_FILE "${_run}/cfl/run.log" ERROR_FILE "${_run}/cfl/error.log")
file(READ "${_run}/cfl/run.log" out)
file(READ "${_run}/cfl/error.log" err)
if(rc EQUAL 0 OR out MATCHES "\\[Step " OR NOT err MATCHES "accepted_step=0 accepted_time=0" OR NOT err MATCHES "half a cell")
    message(FATAL_ERROR "Runtime CFL rejection lost accepted transaction: ${rc}, ${out}${err}")
endif()
file(WRITE "${_run}/verification.json" "{\"status\":\"PASS\",\"smoke_cases\":${_valid},\"accepted_steps\":48,\"regrids\":36,\"invalid_inputs\":${_invalid},\"contextual_ranges\":${_context},\"missing_files\":${_missing},\"extra_input_rejections\":1,\"runtime_rejections\":1}\n")
find_program(POWERSHELL_EXECUTABLE NAMES powershell pwsh REQUIRED)
execute_process(COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
    -File "${CMAKE_CURRENT_LIST_DIR}/interface_first_jet_output_contract.ps1" -Program "${LAMP_EXE}" -RunRoot "${_run}/output_contract"
    RESULT_VARIABLE output_rc OUTPUT_FILE "${_run}/output_contract.log" ERROR_FILE "${_run}/output_contract-error.log")
if(NOT output_rc EQUAL 0)
    message(FATAL_ERROR "FirstJet output contract failed; inspect ${_run}/output_contract-error.log")
endif()
message(STATUS "Five-file FirstJet App contract PASS: ${_valid} smokes; ${_invalid} rejections; cross-file source context, output and rejection checked")
