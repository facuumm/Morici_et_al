% Interpolate if required (2D)
function zint = Interpolate2(x,y,z,valid,mode,maxDistance)

if strcmp(mode,'discard'),
	% In discard mode, do nothing
	zint = z;
else
	% In interpolation mode, interpolate missing points (where time < minTime) using other points
	d = DistanceTransform(valid);
    xx = repmat(x',1,length(y));
    yy = repmat(y,length(x),1);
	if exist('scatteredInterpolant') == 2,
		F = scatteredInterpolant(xx(d==0),yy(d==0),z(d==0));
        zint = F(xx,yy);
    else
        if any(imag(z(:))),
            Freal = TriScatteredInterp(xx(d==0),yy(d==0),real(z(d==0)));
            zintReal = Freal(xx,yy);
            Fimaginary = TriScatteredInterp(xx(d==0),yy(d==0),imag(z(d==0)));
            zintImaginary = Fimaginary(xx,yy);
            zint = zintReal + 1i.*zintImaginary;
        else
            F = TriScatteredInterp(xx(d==0),yy(d==0),z(d==0));
            zint = F(xx,yy);
        end
    end
	% (do not interpolate missing points too distant from valid points)
	zint(d>maxDistance) = z(d>maxDistance);
	zint(isnan(zint)) = z(isnan(zint));
end