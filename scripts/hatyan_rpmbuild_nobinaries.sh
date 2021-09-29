# usage: in PuTTY from current folder: './hatyan_rpmbuild_nobinaries.sh'
# do not forget to chmod +x ./hatyan_rpmbuild.sh
#
# if it doesn't work (permission issue), try:
#    dos2unix hatyan_rpmbuild.sh
#    ./hatyan_rpmbuild.sh
#
# rpmbuild requires (sudo yum -y install): centos-release-scl-rh, rh-python36-python,  rh-python36-python-virtualenv, rpm-build
#!/bin/bash
set -e

#for h6
HATYANROOTFOLDER=~/hatyan_github
RPMTOPDIR=/u/veenstra/rpmbuild #(cannot contain ~ character) is default location on h6 so not per se necessary here

HATYANENVDIR=~/hatyan_env
HATYANEXEC=~/hatyan_fromhome.sh

#delete resulting directories first to start clean
rm -rf ${RPMTOPDIR}
rm -rf ${HATYANENVDIR}

rpmbuild -v -bi ${HATYANROOTFOLDER}/scripts/hatyan_python-latest.spec --define "_topdir ${RPMTOPDIR}" --define "HATYANROOTFOLDER ${HATYANROOTFOLDER}"

cp -r ${RPMTOPDIR}/BUILDROOT/hatyan_python-*/opt/hatyan_python/hatyan_env $HATYANENVDIR
cp ${HATYANROOTFOLDER}/scripts/hatyan.sh $HATYANEXEC
chmod +x $HATYANEXEC
#replace env location in activate script and hatyanexec (double quote to expand variable)
sed -i "s#/opt/hatyan_python/hatyan_env#${HATYANENVDIR}#g" $HATYANENVDIR/bin/activate
sed -i "s#/opt/hatyan_python/hatyan_env#${HATYANENVDIR}#g" $HATYANEXEC

echo "If all went well, there is now a hatyan_env copied from the RPM BUILDROOT in $HATYANENVDIR and $HATYANEXEC can be used to execute hatyan"
