 function  [] = KML_googleEarth(Filename,latitude,longitude,altitude)
 % routine to create .KML file
 
 fod = fopen(Filename, 'wb');
 if fod~=-1
     %     fprintf(fod, '<?xml version=\1.0\ encoding=\UTF-8\?>\n');
     fprintf(fod, '<?xml version="1.0" encoding="UTF-8"?>\n');
     %     fprintf(fod, '<kml xmlns=\http://earth.google.com/kml/2.1\>\n');
     fprintf(fod, '<kml xmlns="http://earth.google.com/kml/2.1">\n');
     fprintf(fod, '<Document>\n');
     %     fprintf(fod, '<Style id=\defaultPlacemark\>\n');
     fprintf(fod, '<Style id="defaultPlacemark">\n');
     fprintf(fod, '\t<IconStyle>\n');
     fprintf(fod, '\t\t<Icon>\n');
     fprintf(fod, '\t\t\t<href>http://maps.google.com/mapfiles/kml/pal4/icon26.png</href>\n');
     fprintf(fod, '\t\t</Icon>\n');
     fprintf(fod, '\t</IconStyle>\n');
     fprintf(fod, '</Style>\n\n');
     %    fprintf(fod, '<Placemark>');
     %    fprintf(fod, '\t<name>Antenna</name>\n');
     %    fprintf(fod, '\t<LookAt>\n');
     %    fprintf(fod, '\t\t<longitude>7.658976358255583</longitude>\n');
     %    fprintf(fod, '\t\t<latitude>45.0652823489628</latitude>\n');
     %    fprintf(fod, '\t\t<altitude>312.0</altitude>\n');
     %    fprintf(fod, '\t</LookAt>\n');
     %    fprintf(fod, '\t<styleUrl>#msn_ylw-pushpin</styleUrl>\n');
     %    fprintf(fod, '\t<Point>\n');
     %    fprintf(fod, '\t\t<coordinates>7.658976358255583,45.0652823489628,312.0</coordinates>\n');
     %    fprintf(fod, '\t</Point>');
     %    fprintf(fod, '</Placemark>\n');
     
     % Reads all file;
     
     %     [values,count] =fread(fid,30000,'float64');
     
     % trova la latitudine longitudine e altezza

     pos_lat=latitude;
     pos_long=longitude;
     pos_height=altitude;
     
     [val] = find(pos_long~= 0);
     
     
     for n=1:length(val)
         fprintf(fod, '<Placemark><styleUrl>#defaultPlacemark</styleUrl>\n');
         fprintf(fod, '<name>%d</name>\n',n);

         fprintf(fod,'<Point>\n\t\t<coordinates>\n');
         %fprintf(fod,'<name>Simple placemark</name>');
         
         fprintf(fod,'\t\t\t%.20f, %.20f, %.20f\n ',pos_long(n),pos_lat(n),pos_height(n));
         fprintf(fod, '\t\t</coordinates>\n\t</Point>\n</Placemark>\n');
     end
     
     fprintf(fod, '\n\n</Document>\n');
     fprintf(fod, '</kml>\n');
     %    fclose(fid);
     fclose(fod);
     
     %     system('PAUSE');
     % return
     %     return 0;
     % }
 else
    errordlg('Input file cannot be opened!','ERROR!');
    error('Err:FileNotOpen','%s file cannot be opened/created.', Filename);
 end
 
 end