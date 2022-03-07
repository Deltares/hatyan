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
#in case of python 3.7/3.8/3.9 instead of 3.6.12 (might be not needed anymore if glibc>=2.14 at RWS)
sed -i "s/python -m pip install pyqt5==5.7.1/python -m pip install pyqt5>=5.7.1/g" hatyan_github/scripts/hatyan_python-latest_python3.spec #unfix old pyqt5 version
sed -i "s# -r %{_topdir}/BUILD/hatyan_github/requirements_dev.txt##g" hatyan_github/scripts/hatyan_python-latest_python3.spec #unfix old dependency versions, since they are not all available in python >=3.8

# setup conda env with specific python version to use in rpmbuild, rpmbuild from spec, deactivate conda env
module load anaconda3 # use ``conda init`` after this command and restart shell, this activates base enviroment and makes it possible to activate venv (report issue to dsc)
conda create -n hatyan_setup_venv python=3.7 -y # 3.6.12 3.7 3.8 3.9 3.10
conda activate hatyan_setup_venv
rpmbuild -v -bb hatyan_github/scripts/hatyan_python-latest_python3.spec --define "VERSIONTAG ${versiontag}"
conda deactivate
module unload anaconda3
echo "RPM was created: $(find $HOME/rpmbuild/RPMS/x86_64/*.rpm | head -n 1)"
