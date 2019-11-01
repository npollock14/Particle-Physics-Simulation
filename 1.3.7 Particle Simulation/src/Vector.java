public class Vector {
	public Direction direction;
	public double magnitude;
	public String tag = "";

	public Vector(Direction dir, double mag) {
		this.direction = dir;
		this.magnitude = mag;
	}

	public Vector(Direction dir, double mag, String tag) {
		this.direction = dir;
		this.magnitude = mag;
		this.tag = tag;
	}
}