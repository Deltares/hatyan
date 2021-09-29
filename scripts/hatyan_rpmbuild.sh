# usage: in PuTTY from current folder: './hatyan_rpmbuild.sh'
# do not forget to chmod +x ./hatyan_rpmbuild.sh
#
# if it doesn't work (permission issue), try:
#    dos2unix hatyan_rpmbuild.sh
#    ./hatyan_rpmbuild.sh
#
# rpmbuild requires (sudo yum -y install): centos-release-scl-rh, rh-python36-python,  rh-python36-python-virtualenv, rpm-build
#!/bin/bash
set -e

#for teamcity
#HATYANROOTFOLDER=%system.teamcity.build.checkoutDir% #workingDir/[definedcheckoutdir]
#RPMTOPDIR=%system.teamcity.build.workingDir%/rpmbuild #(cannot contain ~ character)

#for h6
HATYANROOTFOLDER=~/hatyan_github
RPMTOPDIR=/u/veenstra/rpmbuild #(cannot contain ~ character) is default location on h6 so not per se necessary here

rm -rf ${RPMTOPDIR} #to be sure all RPM's are removed, so quering version number only results in one number
rpmbuild -v -bb ${HATYANROOTFOLDER}/scripts/hatyan_python-latest.spec --define "_topdir ${RPMTOPDIR}" --define "HATYANROOTFOLDER ${HATYANROOTFOLDER}"

#copy specfile and RPM to build folder
HATYANVERSION=$(rpm -q -p ${RPMTOPDIR}/RPMS/x86_64/*.rpm --queryformat '%{VERSION}')
mkdir -p ${HATYANROOTFOLDER}/build
cp ${HATYANROOTFOLDER}/scripts/hatyan_python-latest.spec ${HATYANROOTFOLDER}/build/hatyan_python-${HATYANVERSION}.spec
cp ${RPMTOPDIR}/RPMS/x86_64/*.rpm ${HATYANROOTFOLDER}/build
