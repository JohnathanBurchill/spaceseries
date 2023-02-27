# spaceseries: Docker/Dockerfile

# Copyright (C) 2023  Johnathan K Burchill

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

FROM johnathanburchill/quasidipole:latest AS quasidipole

FROM johnathanburchill/utils:cdev AS dev
USER root
COPY --from=quasidipole /quasidipole.h /usr/local/include/.
COPY --from=quasidipole /libquasidipole.a /usr/local/lib/.
RUN ldconfig
RUN cd /home/science
USER science
WORKDIR /home/science/src
VOLUME ["/home/science/src"]
ENTRYPOINT ["/bin/bash"]



