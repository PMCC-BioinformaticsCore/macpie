document.addEventListener("DOMContentLoaded", function() {
  const toc = document.getElementById("toc");
  if (!toc) return;
  toc.querySelector("h2").addEventListener("click", () => {
    toc.classList.toggle("collapsed");
  });
});

