#!/usr/bin/sh

rm -r b4rpipe.egg-info
rm -r build
rm -r dist

python setup.py bdist_wheel
twine upload --repository pypi dist/*
