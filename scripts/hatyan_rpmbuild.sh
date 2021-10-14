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

versiontag=put_versiontag_here #the versiontag is also internally stored in the specfile, this should be aligned with this one. Possible are main, branches, tags like v2.3.0

#define and delete resulting directories first to start clean
#RPMTOPDIR=%system.teamcity.build.workingDir%/rpmbuild #(cannot contain ~ character) # for teamcity
RPMTOPDIR=/u/veenstra/rpmbuild # for h6: (cannot contain ~ character) is default location on h6 so not per se necessary here
rm -rf ${RPMTOPDIR} #to be sure all RPM's are removed, so quering version number only results in one number

# download spec from source and rpmbuild from spec
rm -rf hatyan_github
git clone -b ${versiontag} https://github.com/Deltares/hatyan.git hatyan_github 
rpmbuild -v -bb hatyan_github/scripts/hatyan_python-latest_python3.spec --define "VERSIONTAG ${versiontag}"
echo "RPM was created: ${RPMTOPDIR}/RPMS/x86_64/*.rpm"
