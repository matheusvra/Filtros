function VectorOut = integral_trap(Vector,TimeSampling)

VectorOut = zeros(1,length(Vector));
VectorOut(1) = Vector(1);
N = length(Vector);
%Loop that perform the integration
for i=2:N                         
	%calculate the integrated vector by trapezoidal method
    VectorOut(i) = VectorOut(i-1) + (Vector(i)+Vector(i-1));
end

%Remove the bias again
%Media = mean(VectorOut)
%VectorOut = VectorOut - Media;
%VectorOut = detrend(VectorOut);

%Adjust the gain
VectorOut = 0.5 * TimeSampling * (VectorOut);
end