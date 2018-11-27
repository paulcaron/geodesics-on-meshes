public class Pair<S, T> {
	public S s;
	public T t;

	public Pair (S s, T t) {
		this.s = s;
		this.t = t;
	}

	public S first() {
		return this.s;
	}

	public T second() {
		return this.t;
	}
}
