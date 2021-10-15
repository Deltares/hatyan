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

#delete resulting directories first to start clean
rm -rf $HOME/rpmbuild #to be sure all RPM's are removed, so quering version number only results in one number. For teamcity it would be %system.teamcity.build.workingDir%/rpmbuild

# download spec from source and rpmbuild from spec
rm -rf hatyan_github
git clone -b ${versiontag} https://github.com/Deltares/hatyan.git hatyan_github 
module load anaconda3
conda create -n hatyan_setup_venv python=3.6.12 -y 
conda activate hatyan_setup_venv
rpmbuild -v -bb hatyan_github/scripts/hatyan_python-latest_python3.spec --define "VERSIONTAG ${versiontag}"
conda deactivate
module unload anaconda3
echo "RPM was created: $(find $HOME/rpmbuild/RPMS/x86_64/*.rpm | head -n 1)"
