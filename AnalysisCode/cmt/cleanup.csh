# echo "cleanup AnalysisCode v0 in /afs/ihep.ac.cn/users/h/huyuxiang/junofs/trunk/offline/Simulation/DetSimV2"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J20v2r0-branch/ExternalLibs/CMT/v1r26
endif
source ${CMTROOT}/mgr/setup.csh
set cmtAnalysisCodetempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtAnalysisCodetempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=AnalysisCode -version=v0 -path=/afs/ihep.ac.cn/users/h/huyuxiang/junofs/trunk/offline/Simulation/DetSimV2  $* >${cmtAnalysisCodetempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=AnalysisCode -version=v0 -path=/afs/ihep.ac.cn/users/h/huyuxiang/junofs/trunk/offline/Simulation/DetSimV2  $* >${cmtAnalysisCodetempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtAnalysisCodetempfile}
  unset cmtAnalysisCodetempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtAnalysisCodetempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtAnalysisCodetempfile}
unset cmtAnalysisCodetempfile
exit $cmtcleanupstatus

