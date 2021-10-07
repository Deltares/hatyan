# usage: in PuTTY from current folder: './hatyan_rpmbuild_nobinaries.sh'
# do not forget to chmod +x ./hatyan_rpmbuild_nobinaries.sh
#
# if it doesn't work (permission issue), try:
#    dos2unix hatyan_rpmbuild.sh
#    ./hatyan_rpmbuild_nobinaries.sh
#
# rpmbuild requires (sudo yum -y install): centos-release-scl-rh, rh-python36-python,  rh-python36-python-virtualenv, rpm-build
#!/bin/bash
set -e

#for h6
RPMTOPDIR=/u/veenstra/rpmbuild #(cannot contain ~ character) is default location on h6 so not per se necessary here
HATYANENVDIR=~/hatyan_env
HATYANEXEC=~/hatyan_fromhome.sh

#delete resulting directories first to start clean
rm -rf ${RPMTOPDIR}
rm -rf ${HATYANENVDIR}

# download spec from source and rpmbuild from spec
versiontag=main #the versiontag is also internally stored in the specfile, this should be aligned with this one. Possible are main, branches, tags like v2.2.68
wget https://github.com/Deltares/hatyan/archive/${versiontag}.zip -O ${versiontag}.zip
rm -rf hatyan-${versiontag} #first delete the destination folder
unzip -p ${versiontag}.zip
rpmbuild -v -bi hatyan-${versiontag}/scripts/hatyan_python-latest_git.spec

cp -r ${RPMTOPDIR}/BUILDROOT/hatyan_python-*/opt/hatyan_python/hatyan_env $HATYANENVDIR
cp hatyan-${versiontag}/scripts/hatyan.sh $HATYANEXEC
chmod +x $HATYANEXEC
#replace env location in activate script and hatyanexec (double quote to expand variable)
sed -i "s#/opt/hatyan_python/hatyan_env#${HATYANENVDIR}#g" $HATYANENVDIR/bin/activate
sed -i "s#/opt/hatyan_python/hatyan_env#${HATYANENVDIR}#g" $HATYANEXEC

echo "If all went well, there is now a hatyan_env copied from the RPM BUILDROOT in $HATYANENVDIR and $HATYANEXEC can be used to execute hatyan"
