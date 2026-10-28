import { redirect } from "next/navigation";

// Legacy route. The old "Impact Explorer" at /results was split into the
// Splicing Impact Explorer (/explore/impact) and the Neoepitope Atlas
// (/explore/neoepitopes). Preserve any ?query= deep links into the impact lens.
export default async function ResultsRedirect({
  searchParams,
}: {
  searchParams: Promise<{ query?: string }>;
}) {
  const params = await searchParams;
  const q = params?.query ? `?query=${encodeURIComponent(params.query)}` : "";
  redirect(`/explore/impact${q}`);
}
