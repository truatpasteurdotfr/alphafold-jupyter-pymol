# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

FROM  ghcr.io/truatpasteurdotfr/alphafold-jupyter:main

# Install conda packages.
RUN conda install -qy \
        msgpack-cxx \
        msgpack-c \
        libgl-devel \
        freetype \
        glew \
        glm \
        libglu \
        libnetcdf \
        libpng \
        libxml2-devel \
        pmw \
        pyqt \
        xorg-libxfixes
RUN	git clone https://github.com/schrodinger/pymol-open-source.git && \
	git clone https://github.com/rcsb/mmtf-cpp.git && \
	mv mmtf-cpp/include/mmtf* pymol-open-source/include/ && \
	cd pymol-open-source && \
	python3 setup.py build install 
RUN	rm -rf pymol-open-source mmtf-cpp
