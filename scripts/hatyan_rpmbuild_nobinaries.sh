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

versiontag=main #the versiontag is also internally stored in the specfile, this should be aligned with this one. Possible are main, branches, tags like v2.2.68

#define and delete resulting directories first to start clean (for h6)
RPMTOPDIR=/u/veenstra/rpmbuild #(cannot contain ~ character) is default location on h6 so not per se necessary here
HATYANENVDIR=~/hatyan_env
HATYANEXEC=~/hatyan_fromhome.sh
rm -rf ${RPMTOPDIR}
rm -rf ${HATYANENVDIR}
rm -f ${HATYANEXEC}

# download spec from source and rpmbuild from spec
rm -rf hatyan_github
git clone -b ${versiontag} https://github.com/Deltares/hatyan.git hatyan_github 
rpmbuild -v -bi hatyan_github/scripts/hatyan_python-latest_git.spec --define "VERSIONTAG main"

cp -r ${RPMTOPDIR}/BUILDROOT/hatyan_python-*/opt/hatyan_python/hatyan_env $HATYANENVDIR
cp hatyan_github/scripts/hatyan.sh $HATYANEXEC
chmod +x $HATYANEXEC
#replace env location in activate script and hatyanexec (double quote to expand variable)
sed -i "s#/opt/hatyan_python/hatyan_env#${HATYANENVDIR}#g" $HATYANENVDIR/bin/activate
sed -i "s#/opt/hatyan_python/hatyan_env#${HATYANENVDIR}#g" $HATYANEXEC

echo "If all went well, there is now a hatyan_env copied from the RPM BUILDROOT in $HATYANENVDIR and $HATYANEXEC can be used to execute hatyan:"
echo "EXAMPLE: ./hatyan_fromhome.sh hatyan-${versiontag}/tests/configfiles/predictie_2019_19Ycomp4Ydia_VLISSGN_interactive.py"
