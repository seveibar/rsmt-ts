import test from "ava"
import { rfst } from "../index"

test("rectilinear steiner minimal tree", async (t) => {
  const nodes: [number, number][] = [
    [0, 0],
    [1, 2],
    [4, 1],
  ]
  const solution = rfst(nodes)
  t.snapshot(solution)
})
