# echo "setup AnalysisCode v0 in /junofs/users/huyuxiang/juno_centos7/offline/Simulation/DetSimV2"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J20v2r0-branch/ExternalLibs/CMT/v1r26
endif
source ${CMTROOT}/mgr/setup.csh
set cmtAnalysisCodetempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtAnalysisCodetempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=AnalysisCode -version=v0 -path=/junofs/users/huyuxiang/juno_centos7/offline/Simulation/DetSimV2  -no_cleanup $* >${cmtAnalysisCodetempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=AnalysisCode -version=v0 -path=/junofs/users/huyuxiang/juno_centos7/offline/Simulation/DetSimV2  -no_cleanup $* >${cmtAnalysisCodetempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtAnalysisCodetempfile}
  unset cmtAnalysisCodetempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtAnalysisCodetempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtAnalysisCodetempfile}
unset cmtAnalysisCodetempfile
exit $cmtsetupstatus

