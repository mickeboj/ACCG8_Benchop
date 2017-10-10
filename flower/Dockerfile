FROM lordnelsson/accbop_base
ADD /BENCHOP bench/
ADD /cproj cproj/
ADD /flower/config.py cproj/
RUN pip install flower
EXPOSE 5555
CMD flower -A cproj --port=5555
