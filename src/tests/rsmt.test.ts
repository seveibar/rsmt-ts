import test from "ava"
import rsmt from "../index"

test("rectilinear steiner minimal tree", async (t) => {
  const nodes = [
    [0, 0],
    [1, 2],
    [4, 1],
  ]
  const solution = await rsmt(nodes)
  t.deepEqual(solution.edges, [
    [
      [0, 0],
      [0, 1],
    ],
    [
      [1, 1],
      [1, 2],
    ],
    [
      [1, 1],
      [4, 1],
    ],
    [
      [1, 1],
      [0, 1],
    ],
  ])
})
