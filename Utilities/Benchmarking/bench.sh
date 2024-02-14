#!/usr/bin/env bash

for file in build/bin/*
do
	output=$(./${file})
	keygen_time=0
	verification_time=0
	signature_time=0
	while IFS= read -r line
	do
		if [[ "${line}" == *"Key generation kCycles"* ]]; then
			times=$(echo ${line} | sed 's/[^0-9]*//')
			times=(${times//,/ })
			kg_time=${times[0]}
		fi
		if [[ "${line}" == *"Signature kCycles"* ]]; then
			times=$(echo ${line} | sed 's/[^0-9]*//')
			times=(${times//,/ })
			signature_time=${times[0]}
			# echo ${signature_time}
		fi
		if [[ "${line}" == *"Verification kCycles"* ]]; then
			times=$(echo ${line} | sed 's/[^0-9]*//')
			times=(${times//,/ })
			verification_time=${times[0]}
			# echo ${verification_time}
		fi
	done <<< "${output}"

    echo "${file} | ${kg_time} | ${signature_time} | ${verification_time}"
done
