function kod=tablicaKodera( kod, tr )
if( 0==isempty( tr.value ) )
    kod(length(kod)+1).value = tr.value;
    kod(end).bits = [];
    return;
end
kodleft = tablicaKodera( kod, tr.lewy );
kodright = tablicaKodera( kod, tr.prawy );
for n=1:length(kodleft)
    kodleft(n).bits = ['1', kodleft(n).bits];
end
for n=1:length(kodright)
    kodright(n).bits = ['0', kodright(n).bits];
end
kod = [kod, kodleft, kodright];