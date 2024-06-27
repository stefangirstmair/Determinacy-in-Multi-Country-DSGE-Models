function ax=loadTile(tc, name, ID)
 p=openfig(name,'new','invisible');
 ax=findobj(p,'Type','axes');
 ax.Parent=tc;
 ax.Layout.Tile=ID;
 close(p)
end