#!/bin/csh

# run this script after you compile the translate program

# directory containing working copies of the tahoe modules
setenv TAHOE /opt/sandia

# username on Tahoe server
set TAHOE_USER = sawimme

# host name - used only for subject line of test report
set HOST = arisun.nrl.navy.mil

# paths to commands
set ECHO = '/usr/bin/echo'
set MAKE = '/usr/ccs/bin/make'

######### hopefully no changes will be required below #########

set BENCHDIR = ${TAHOE}/contrib/translate/benchmarks

# clean old benchmark results
${ECHO} "removing old benchmark tests"
cd ${BENCHDIR}; ${MAKE} clean

# run benchmarks
${ECHO} "running benchmarks"
cd ${BENCHDIR}; ${MAKE} trans

# run comparator
${ECHO} "running comparator"
    if(-e ${BENCHDIR}/compare.summary) then
	rm ${BENCHDIR}/compare.summary
    endif
    if(-e ${BENCHDIR}/compare.console) then
	rm ${BENCHDIR}/compare.console
    endif
cd ${BENCHDIR}; ${TAHOE}/benchmark/comparator/compare -f compare.batch >&! compare.console

# create summary
${ECHO} "creating summary"
    if(! -e ${BENCHDIR}/compare.summary) then
#       test/s failed - did not complete
	${ECHO} "BENCHMARKs failed"
		${ECHO} "Subject: TAHOE: regression FAIL ${HOST}\nTAHOE regression FAIL:" `date` > ${BENCHDIR}/tmp
		cat ${BENCHDIR}/tmp ${BENCHDIR}/compare.console > ${BENCHDIR}/tmp_mail
    else if(`grep -c "FAIL" ${BENCHDIR}/compare.summary` ) then
#       test/s failed
	${ECHO} "BENCHMARKs failed"
		${ECHO} "Subject: TAHOE: regression FAIL ${HOST}\nTAHOE regression FAIL:" `date` > ${BENCHDIR}/tmp
		cat ${BENCHDIR}/tmp ${BENCHDIR}/compare.summary > ${BENCHDIR}/tmp_mail
    else
	${ECHO} "BENCHMARKs passed"
		${ECHO} "Subject: TAHOE: regression PASS ${HOST}\nTAHOE regression PASS:" `date` > ${BENCHDIR}/tmp
		cat ${BENCHDIR}/tmp ${BENCHDIR}/compare.summary > ${BENCHDIR}/tmp_mail
	endif
    else
#       make failed
		${ECHO} "Subject: TAHOE: build FAIL ${HOST}\nTAHOE build FAIL:" `date` > ${BENCHDIR}/tmp
		cat ${BENCHDIR}/tmp ${TAHOE}/make.stderr > ${BENCHDIR}/tmp_mail
    endif
