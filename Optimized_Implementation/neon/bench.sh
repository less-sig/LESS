#!/usr/bin/env sh

calc(){ awk "BEGIN { print "$*" }"; }

for file in cmake-build-release/*
do
    if grep -q benchmark "${file}"; then
	    output=$(sudo ./${file} 2>&1)
	    keygen_time=0
	    verification_time=0
	    signature_time=0
	    while IFS= read -r line
	    do
	    	if [[ "${line}" == *"Key generation kCycles"* ]]; then
	    		times=$(echo ${line} | sed 's/[^0-9]*//')
	    		times=(${times//,/ })
                kg_time=$(calc ${times[0]}/1000)
	    		#kg_time=${times[0]}
	    	fi
	    	if [[ "${line}" == *"Signature kCycles"* ]]; then
	    		times=$(echo ${line} | sed 's/[^0-9]*//')
	    		times=(${times//,/ })
                signature_time=$(calc ${times[0]}/1000)
	    		#signature_time=${times[0]}
	    	fi
	    	if [[ "${line}" == *"Verification kCycles"* ]]; then
	    		times=$(echo ${line} | sed 's/[^0-9]*//')
	    		times=(${times//,/ })
                verification_time=$(calc ${times[0]}/1000)
	    		#verification_time=${times[0]}
	    	fi
	    done <<< "${output}"

        echo "${file} | ${kg_time} | ${signature_time} | ${verification_time}"
    fi
done
