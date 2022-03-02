elevation_ref = load("x_elevation_ref.mat");
travel_ref = load("x_travel_ref.mat");

elevation_ref.ans(4)
plot(elevation_ref.ans(1,:), elevation_ref.ans(2,:))