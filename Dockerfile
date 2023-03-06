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
ENV PATH="/opt/conda/bin:$PATH"
RUN apt-get update && apt-get -y install \
	gcc            \
	libglew-dev    \
	libglm-dev     \
	libnetcdf-dev  \
	libmsgpack-dev \
	libxi6 \
	libxinerama1 \
	xkb-data \
	libxkbcommon0
RUN	git clone https://github.com/schrodinger/pymol-open-source.git && \
	git clone https://github.com/rcsb/mmtf-cpp.git && \
	mv mmtf-cpp/include/mmtf* pymol-open-source/include/ && \
	cd pymol-open-source && \
	python3 setup.py build install 
RUN	rm -rf pymol-open-source mmtf-cpp
# pymol runtime
RUN	apt-get update &&  \
	DEBIAN_FRONTEND=noninteractive apt-get -y install \
	libdbus-1-3 \ 
	libfontconfig1 \ 
	libgl1-mesa-glx \ 
	libglib2.0-0 \ 
	libice6 \ 
	libsm6 \ 
	libxau6 \ 
	libxcb-icccm4 \  
	libxcb-image0 \ 
	libxcb-keysyms1 \ 
	libxcb-randr0 \ 
	libxcb-render-util0 \ 
	libxcb-shape0 \ 
	libxcb-xfixes0 \ 
	libxcb-xinerama0 \ 
	libxcb-xinput0 \ 
	libxi6 \ 
	libxkbcommon-x11-0 \ 
	libxrender1 \ 
	mesa-utils 

